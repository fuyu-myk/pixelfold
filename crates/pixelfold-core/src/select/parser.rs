//! The nom parser: a query string to a [`Selection`] tree.
//!
//! Recursive descent over the grammar in the module docs. Precedence falls out
//! of the call chain: `or` calls `and` calls `not` calls `primary`.

use nom::branch::alt;
use nom::bytes::complete::{tag, take_while1};
use nom::character::complete::{char, multispace0, u32 as u32_of, u64 as u64_of};
use nom::combinator::{all_consuming, opt, value, verify};
use nom::multi::many1;
use nom::number::complete::float;
use nom::sequence::{delimited, preceded};
use nom::{IResult, Parser};

use super::{Cmp, ResiRange, SelectError, Selection};
use crate::structure::SecondaryStructure;

/// The deepest a selection may nest. A query is recursive descent, so its
/// nesting is bounded to keep a crafted `((((...))))` or `not not not ...` from
/// overflowing the stack; this also bounds the evaluator, which recurses over
/// the same tree. Any real selection is a handful of levels deep.
const MAX_DEPTH: u32 = 128;

/// Parse a complete query, rejecting any trailing input.
pub fn parse(query: &str) -> Result<Selection, SelectError> {
    let top = |input| selection(input, 0);
    match all_consuming(delimited(multispace0, top, multispace0)).parse(query) {
        Ok((_, selection)) => Ok(selection),
        Err(_) => Err(SelectError {
            message: format!("could not parse selection: {query:?}"),
        }),
    }
}

/// A hard parse failure: `Failure` stops parsing rather than backtracking, so a
/// depth cap or an out-of-range number is reported at once, not retried down
/// another branch.
fn reject(input: &str) -> nom::Err<nom::error::Error<&str>> {
    nom::Err::Failure(nom::error::Error::new(input, nom::error::ErrorKind::Verify))
}

/// The reserved words. A bare name in a list stops at any of them, so an
/// operator or another predicate after a list is never swallowed into it.
fn is_keyword(word: &str) -> bool {
    matches!(
        word.to_ascii_lowercase().as_str(),
        "and"
            | "or"
            | "not"
            | "of"
            | "byres"
            | "within"
            | "all"
            | "none"
            | "chain"
            | "resi"
            | "resn"
            | "name"
            | "element"
            | "model"
            | "ss"
            | "bfactor"
            | "plddt"
            | "protein"
            | "nucleic"
            | "ligand"
            | "water"
            | "metal"
            | "hetatm"
    )
}

/// A bare token: letters, digits, and the primes and stars of nucleic atom
/// names. Ranges and comparisons are parsed elsewhere, so no punctuation here.
fn token(input: &str) -> IResult<&str, &str> {
    preceded(
        multispace0,
        take_while1(|c: char| c.is_alphanumeric() || c == '\'' || c == '*'),
    )
    .parse(input)
}

/// Match one whole reserved word, so `or` does not match the start of `organic`
/// and `ss` does not match `ssbond`.
fn keyword<'a>(word: &'static str) -> impl FnMut(&'a str) -> IResult<&'a str, ()> {
    move |input| {
        let (rest, matched) = token(input)?;
        if matched.eq_ignore_ascii_case(word) {
            Ok((rest, ()))
        } else {
            Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Tag,
            )))
        }
    }
}

/// One or more bare names, each not a reserved word.
fn name_list(input: &str) -> IResult<&str, Vec<String>> {
    many1(verify(token, |t: &str| !is_keyword(t)))
        .parse(input)
        .map(|(rest, names)| (rest, names.into_iter().map(String::from).collect()))
}

fn selection(input: &str, depth: u32) -> IResult<&str, Selection> {
    or_expr(input, depth)
}

fn or_expr(input: &str, depth: u32) -> IResult<&str, Selection> {
    let (mut rest, mut acc) = and_expr(input, depth)?;
    while let Ok((next, _)) = keyword("or")(rest) {
        let (next, rhs) = and_expr(next, depth)?;
        acc = Selection::Or(Box::new(acc), Box::new(rhs));
        rest = next;
    }

    Ok((rest, acc))
}

fn and_expr(input: &str, depth: u32) -> IResult<&str, Selection> {
    let (mut rest, mut acc) = not_expr(input, depth)?;
    while let Ok((next, _)) = keyword("and")(rest) {
        let (next, rhs) = not_expr(next, depth)?;
        acc = Selection::And(Box::new(acc), Box::new(rhs));
        rest = next;
    }

    Ok((rest, acc))
}

fn not_expr(input: &str, depth: u32) -> IResult<&str, Selection> {
    if depth >= MAX_DEPTH {
        return Err(reject(input));
    }

    if let Ok((rest, _)) = keyword("not")(input) {
        let (rest, inner) = not_expr(rest, depth + 1)?;
        return Ok((rest, Selection::Not(Box::new(inner))));
    }

    primary(input, depth)
}

fn primary(input: &str, depth: u32) -> IResult<&str, Selection> {
    if depth >= MAX_DEPTH {
        return Err(reject(input));
    }

    alt((
        |i| parenthesised(i, depth),
        |i| byres(i, depth),
        |i| within(i, depth),
        predicate,
    ))
    .parse(input)
}

fn parenthesised(input: &str, depth: u32) -> IResult<&str, Selection> {
    delimited(
        preceded(multispace0, char('(')),
        |i| selection(i, depth + 1),
        preceded(multispace0, char(')')),
    )
    .parse(input)
}

fn byres(input: &str, depth: u32) -> IResult<&str, Selection> {
    let (rest, _) = keyword("byres")(input)?;
    let (rest, inner) = primary(rest, depth + 1)?;
    Ok((rest, Selection::ByRes(Box::new(inner))))
}

fn within(input: &str, depth: u32) -> IResult<&str, Selection> {
    let (rest, _) = keyword("within")(input)?;
    let (rest, distance) = preceded(multispace0, float).parse(rest)?;
    if !distance.is_finite() || distance < 0.0 {
        return Err(reject(rest)); // a distance must be a finite, non-negative number
    }

    let (rest, _) = keyword("of")(rest)?;
    let (rest, inner) = primary(rest, depth + 1)?;

    Ok((rest, Selection::Within(distance, Box::new(inner))))
}

fn predicate(input: &str) -> IResult<&str, Selection> {
    alt((nullary, listed, numeric)).parse(input)
}

/// The predicates that take no argument.
fn nullary(input: &str) -> IResult<&str, Selection> {
    alt((
        value(Selection::All, keyword("all")),
        value(Selection::None, keyword("none")),
        value(Selection::Protein, keyword("protein")),
        value(Selection::Nucleic, keyword("nucleic")),
        value(Selection::Ligand, keyword("ligand")),
        value(Selection::Water, keyword("water")),
        value(Selection::Metal, keyword("metal")),
        value(Selection::Hetatm, keyword("hetatm")),
    ))
    .parse(input)
}

/// The predicates that take a list of names.
fn listed(input: &str) -> IResult<&str, Selection> {
    alt((
        named("chain", Selection::Chain),
        named("resn", Selection::Resn),
        named("name", Selection::Name),
        named("element", Selection::Element),
    ))
    .parse(input)
}

fn named<'a>(
    word: &'static str,
    build: fn(Vec<String>) -> Selection,
) -> impl FnMut(&'a str) -> IResult<&'a str, Selection> {
    move |input| {
        let (rest, _) = keyword(word)(input)?;
        let (rest, names) = name_list(rest)?;
        Ok((rest, build(names)))
    }
}

/// The predicates that take numbers or a comparison.
fn numeric(input: &str) -> IResult<&str, Selection> {
    alt((resi, model, secondary, bfactor, plddt)).parse(input)
}

fn resi(input: &str) -> IResult<&str, Selection> {
    let (rest, _) = keyword("resi")(input)?;
    let (rest, ranges) = many1(preceded(multispace0, resi_range)).parse(rest)?;
    Ok((rest, Selection::Resi(ranges)))
}

fn resi_range(input: &str) -> IResult<&str, ResiRange> {
    let (rest, lo) = u32_of(input)?;
    let (rest, hi) = opt(preceded(char('-'), u32_of)).parse(rest)?;
    Ok((
        rest,
        ResiRange {
            lo,
            hi: hi.unwrap_or(lo),
        },
    ))
}

fn model(input: &str) -> IResult<&str, Selection> {
    let (rest, _) = keyword("model")(input)?;
    let (rest, number) = preceded(multispace0, u64_of).parse(rest)?;
    Ok((rest, Selection::Model(number as usize)))
}

fn secondary(input: &str) -> IResult<&str, Selection> {
    let (rest, _) = keyword("ss")(input)?;
    let (rest, letter) = token(rest)?;
    let ss = match letter.to_ascii_uppercase().as_str() {
        "H" => SecondaryStructure::Helix,
        "E" => SecondaryStructure::Sheet,
        "T" => SecondaryStructure::Turn,
        "C" => SecondaryStructure::Coil,
        _ => {
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Tag,
            )));
        }
    };

    Ok((rest, Selection::Ss(ss)))
}

fn bfactor(input: &str) -> IResult<&str, Selection> {
    let (rest, _) = keyword("bfactor")(input)?;
    let (rest, (cmp, bound)) = comparison(rest)?;
    Ok((rest, Selection::BFactor(cmp, bound)))
}

fn plddt(input: &str) -> IResult<&str, Selection> {
    let (rest, _) = keyword("plddt")(input)?;
    let (rest, (cmp, bound)) = comparison(rest)?;
    Ok((rest, Selection::Plddt(cmp, bound)))
}

fn comparison(input: &str) -> IResult<&str, (Cmp, f32)> {
    let (rest, cmp) = preceded(
        multispace0,
        alt((
            value(Cmp::Le, tag("<=")),
            value(Cmp::Ge, tag(">=")),
            value(Cmp::Eq, tag("==")),
            value(Cmp::Eq, tag("=")),
            value(Cmp::Lt, tag("<")),
            value(Cmp::Gt, tag(">")),
        )),
    )
    .parse(input)?;
    let (rest, bound) = preceded(multispace0, float).parse(rest)?;

    Ok((rest, (cmp, bound)))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parsed(query: &str) -> Selection {
        parse(query).unwrap()
    }

    #[test]
    fn simple_predicates_parse() {
        assert_eq!(parsed("all"), Selection::All);
        assert_eq!(parsed("none"), Selection::None);
        assert_eq!(parsed("water"), Selection::Water);
        assert_eq!(
            parsed("chain A B"),
            Selection::Chain(vec!["A".into(), "B".into()])
        );
        assert_eq!(parsed("resn ALA"), Selection::Resn(vec!["ALA".into()]));
        assert_eq!(parsed("model 2"), Selection::Model(2));
        assert_eq!(parsed("ss H"), Selection::Ss(SecondaryStructure::Helix));
    }

    #[test]
    fn keywords_are_case_insensitive_and_whole_words() {
        assert_eq!(parsed("ALL"), Selection::All);
        assert_eq!(parsed("Chain a"), Selection::Chain(vec!["a".into()]));
        // A name sharing a prefix with a keyword is still a name: "ORN" is not
        // the operator "or".
        assert_eq!(parsed("resn ORN"), Selection::Resn(vec!["ORN".into()]));
        // But a name equal to a reserved word cannot be selected bare.
        assert!(parse("resn all").is_err());
    }

    #[test]
    fn residue_ranges_and_singletons() {
        assert_eq!(
            parsed("resi 1-10 20 30-40"),
            Selection::Resi(vec![
                ResiRange { lo: 1, hi: 10 },
                ResiRange { lo: 20, hi: 20 },
                ResiRange { lo: 30, hi: 40 },
            ])
        );
    }

    #[test]
    fn comparisons_take_the_longest_operator() {
        assert_eq!(parsed("bfactor <= 30"), Selection::BFactor(Cmp::Le, 30.0));
        assert_eq!(parsed("bfactor < 30"), Selection::BFactor(Cmp::Lt, 30.0));
        assert_eq!(parsed("plddt >= 70"), Selection::Plddt(Cmp::Ge, 70.0));
        assert_eq!(parsed("bfactor = 50"), Selection::BFactor(Cmp::Eq, 50.0));
    }

    #[test]
    fn precedence_is_not_over_and_over_or() {
        // not a and b or c  ==  ((not a) and b) or c
        let parsed = parsed("not water and chain A or metal");
        let expected = Selection::Or(
            Box::new(Selection::And(
                Box::new(Selection::Not(Box::new(Selection::Water))),
                Box::new(Selection::Chain(vec!["A".into()])),
            )),
            Box::new(Selection::Metal),
        );
        assert_eq!(parsed, expected);
    }

    #[test]
    fn parentheses_override_precedence() {
        let parsed = parsed("chain A and (resn ALA or resn GLY)");
        let expected = Selection::And(
            Box::new(Selection::Chain(vec!["A".into()])),
            Box::new(Selection::Or(
                Box::new(Selection::Resn(vec!["ALA".into()])),
                Box::new(Selection::Resn(vec!["GLY".into()])),
            )),
        );
        assert_eq!(parsed, expected);
    }

    #[test]
    fn within_and_byres_bind_to_a_primary() {
        assert_eq!(
            parsed("within 5 of resn HEM"),
            Selection::Within(5.0, Box::new(Selection::Resn(vec!["HEM".into()])))
        );
        assert_eq!(
            parsed("byres chain A"),
            Selection::ByRes(Box::new(Selection::Chain(vec!["A".into()])))
        );
        // The operand is a primary, so a compound one needs parentheses.
        assert_eq!(
            parsed("within 4.5 of (metal or water)"),
            Selection::Within(
                4.5,
                Box::new(Selection::Or(
                    Box::new(Selection::Metal),
                    Box::new(Selection::Water)
                ))
            )
        );
    }

    #[test]
    fn a_list_stops_at_the_next_keyword() {
        // resn's list ends at "and", not swallowing the next predicate.
        assert_eq!(
            parsed("resn ALA and chain B"),
            Selection::And(
                Box::new(Selection::Resn(vec!["ALA".into()])),
                Box::new(Selection::Chain(vec!["B".into()])),
            )
        );
    }

    #[test]
    fn malformed_queries_error() {
        assert!(parse("").is_err());
        assert!(parse("chain").is_err()); // a list predicate with no names
        assert!(parse("all extra").is_err()); // trailing input after a complete term
        assert!(parse("(chain A").is_err()); // unbalanced
        assert!(parse("within of metal").is_err()); // missing distance
        assert!(parse("ss X").is_err()); // not a structure class
    }

    #[test]
    fn a_negative_or_infinite_within_distance_is_rejected() {
        assert!(parse("within -5 of metal").is_err());
        assert!(parse("within inf of metal").is_err());
        // A zero distance is allowed: it selects the seed itself.
        assert_eq!(
            parsed("within 0 of metal"),
            Selection::Within(0.0, Box::new(Selection::Metal))
        );
    }

    #[test]
    fn deeply_nested_queries_error_instead_of_overflowing_the_stack() {
        // A crafted query that would recurse past the stack is rejected as a
        // parse error, not a process abort. Well past MAX_DEPTH on every axis.
        let parens = format!("{}all{}", "(".repeat(5000), ")".repeat(5000));
        assert!(parse(&parens).is_err());
        assert!(parse(&format!("{}all", "not ".repeat(5000))).is_err());
        assert!(parse(&format!("{}all", "byres ".repeat(5000))).is_err());
        assert!(parse(&format!("{}all", "within 5 of ".repeat(5000))).is_err());

        // A realistically nested query still parses.
        let nested = format!("{}metal{}", "(".repeat(50), ")".repeat(50));
        assert!(parse(&nested).is_ok());
    }
}
