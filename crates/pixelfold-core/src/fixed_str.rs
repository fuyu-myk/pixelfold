//! Inline fixed-capacity ASCII strings for hot struct fields.
//!
//! Atom names, residue names, and chain ids are short and compared constantly.
//! Storing them as `[u8; N]` keeps `Atom` `Copy`-friendly and allocation-free,
//! and makes them cheap `HashMap` keys. Content longer than `N` bytes (extended
//! chemical-component codes, long mmCIF chain ids) is truncated on construction.

/// A fixed-capacity ASCII string of at most `N` bytes. The value is the bytes up
/// to the first NUL; construction zero-pads and truncates.
///
/// Because construction zero-pads, the derived byte-wise ordering matches
/// lexicographic string ordering (NUL sorts before any printable byte).
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct FixedStr<const N: usize> {
    bytes: [u8; N],
}

impl<const N: usize> FixedStr<N> {
    /// Build from a string slice, keeping the first `N` bytes.
    pub fn new(s: &str) -> Self {
        let mut bytes = [0u8; N];
        let src = s.as_bytes();
        let len = src.len().min(N);
        bytes[..len].copy_from_slice(&src[..len]);

        Self { bytes }
    }

    /// The stored value as a string slice (empty if truncation split a
    /// multi-byte character, which cannot happen for ASCII structure data).
    pub fn as_str(&self) -> &str {
        let end = self.bytes.iter().position(|&b| b == 0).unwrap_or(N);
        std::str::from_utf8(&self.bytes[..end]).unwrap_or("")
    }

    /// True when no bytes are stored.
    pub fn is_empty(&self) -> bool {
        self.bytes[0] == 0
    }
}

impl<const N: usize> std::fmt::Display for FixedStr<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

impl<const N: usize> std::fmt::Debug for FixedStr<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self.as_str())
    }
}

impl<const N: usize> PartialEq<str> for FixedStr<N> {
    fn eq(&self, other: &str) -> bool {
        self.as_str() == other
    }
}

impl<const N: usize> PartialEq<&str> for FixedStr<N> {
    fn eq(&self, other: &&str) -> bool {
        self.as_str() == *other
    }
}

impl<const N: usize> From<&str> for FixedStr<N> {
    fn from(s: &str) -> Self {
        Self::new(s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn round_trips_shorter_and_exact_lengths() {
        assert_eq!(FixedStr::<4>::new("CA").as_str(), "CA");
        assert_eq!(FixedStr::<4>::new("N").as_str(), "N");
        assert_eq!(FixedStr::<4>::new("HD11").as_str(), "HD11"); // full capacity
        assert_eq!(FixedStr::<3>::new("ALA").as_str(), "ALA");
        assert_eq!(FixedStr::<2>::new("A").as_str(), "A");
    }

    #[test]
    fn truncates_beyond_capacity() {
        assert_eq!(FixedStr::<3>::new("TOOLONG").as_str(), "TOO");
        assert_eq!(FixedStr::<2>::new("").as_str(), "");
        assert!(FixedStr::<2>::new("").is_empty());
    }

    #[test]
    fn compares_against_str_and_itself() {
        let ca = FixedStr::<4>::new("CA");
        assert!(ca == "CA");
        assert!(ca != "CB");
        assert_eq!(ca, FixedStr::<4>::new("CA"));
        assert_ne!(ca, FixedStr::<4>::new("N"));
    }

    #[test]
    fn displays_without_padding() {
        assert_eq!(format!("{}", FixedStr::<3>::new("SER")), "SER");
        assert_eq!(format!("{:?}", FixedStr::<3>::new("SER")), "\"SER\"");
    }

    #[test]
    fn usable_as_hash_key() {
        let mut map: HashMap<(FixedStr<2>, u32), u8> = HashMap::new();
        map.insert((FixedStr::new("A"), 10), 1);
        assert_eq!(map.get(&(FixedStr::<2>::new("A"), 10)), Some(&1));
        assert_eq!(map.get(&(FixedStr::<2>::new("B"), 10)), None);
    }
}
