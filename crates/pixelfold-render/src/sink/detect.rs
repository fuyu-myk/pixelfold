//! Terminal graphics capability detection from the environment.
//!
//! Querying the terminal directly is more precise but needs a live TTY in raw
//! mode; the environment covers the common terminals and is what a headless test
//! can exercise. The caller supplies the lookup, so detection stays pure.

/// A terminal image protocol, in order of fidelity.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Protocol {
    Kitty,
    Iterm2,
    HalfBlock,
}

/// Pick the best protocol the environment advertises. `getenv` returns a
/// variable's value, mirroring `std::env::var(..).ok()`.
pub fn detect(getenv: impl Fn(&str) -> Option<String>) -> Protocol {
    let set = |k: &str| getenv(k).is_some();
    let val = |k: &str| getenv(k).unwrap_or_default();
    let term = val("TERM");
    let program = val("TERM_PROGRAM");

    // Kitty's protocol, and the terminals that also implement it.
    if set("KITTY_WINDOW_ID")
        || term.contains("kitty")
        || term.contains("ghostty")
        || program == "ghostty"
        || program == "WezTerm"
    {
        return Protocol::Kitty;
    }

    if program == "iTerm.app" || set("ITERM_SESSION_ID") {
        return Protocol::Iterm2;
    }

    Protocol::HalfBlock
}

/// Whether the session looks remote, where the compressed protocols matter most.
pub fn is_ssh(getenv: impl Fn(&str) -> Option<String>) -> bool {
    getenv("SSH_CONNECTION").is_some() || getenv("SSH_TTY").is_some()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A lookup over a fixed set of pairs.
    fn env<'a>(pairs: &'a [(&'a str, &'a str)]) -> impl Fn(&str) -> Option<String> + 'a {
        move |k| {
            pairs
                .iter()
                .find(|(key, _)| *key == k)
                .map(|(_, v)| v.to_string())
        }
    }

    #[test]
    fn kitty_is_detected_by_window_id_and_by_term() {
        assert_eq!(detect(env(&[("KITTY_WINDOW_ID", "1")])), Protocol::Kitty);
        assert_eq!(detect(env(&[("TERM", "xterm-kitty")])), Protocol::Kitty);
        assert_eq!(detect(env(&[("TERM", "xterm-ghostty")])), Protocol::Kitty);
        assert_eq!(detect(env(&[("TERM_PROGRAM", "WezTerm")])), Protocol::Kitty);
    }

    #[test]
    fn iterm2_is_detected_by_term_program() {
        assert_eq!(
            detect(env(&[("TERM_PROGRAM", "iTerm.app")])),
            Protocol::Iterm2
        );
        assert_eq!(detect(env(&[("ITERM_SESSION_ID", "w0")])), Protocol::Iterm2);
    }

    #[test]
    fn an_unknown_terminal_falls_back_to_half_blocks() {
        assert_eq!(
            detect(env(&[("TERM", "xterm-256color")])),
            Protocol::HalfBlock
        );
        assert_eq!(detect(env(&[])), Protocol::HalfBlock);
    }

    #[test]
    fn ssh_is_detected_from_either_variable() {
        assert!(is_ssh(env(&[("SSH_CONNECTION", "1.2.3.4 5 6.7.8.9 22")])));
        assert!(is_ssh(env(&[("SSH_TTY", "/dev/pts/0")])));
        assert!(!is_ssh(env(&[("TERM", "xterm-kitty")])));
    }
}
