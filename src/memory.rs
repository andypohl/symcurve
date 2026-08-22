//! A memory budget, and the chunk size it implies.

use std::fmt;
use std::str::FromStr;

const KI: u64 = 1024;
const MI: u64 = KI * 1024;
const GI: u64 = MI * 1024;
const TI: u64 = GI * 1024;

/// An upper bound on the memory used for buffering scores.
///
/// This governs how much scored output is held at once, which is the part of the
/// footprint that scales with the genome. It does not cover the sequence of the record
/// being read: without an index the FASTA reader hands over a whole record at a time, so
/// the largest record's sequence is a floor that a budget cannot lower.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MemoryBudget(u64);

/// Bytes held per buffered score. Scores are computed and buffered as `f64`.
const BYTES_PER_SCORE: u64 = 8;

/// Never go below this many scores per chunk; tiny chunks pay more in per-chunk lead-in
/// than they save in memory.
const MIN_CHUNK_SCORES: usize = 4096;

impl MemoryBudget {
    pub fn bytes(&self) -> u64 {
        self.0
    }

    /// How many scores one chunk should produce, given the number of worker threads.
    ///
    /// Each thread holds one chunk's scores while working, so the in-flight cost is
    /// `threads * chunk_scores * 8` bytes. The result is capped at `cap` because past a
    /// point larger chunks only reduce parallelism, and floored so that a very small
    /// budget degrades rather than grinding to a halt.
    pub fn chunk_scores(&self, threads: usize, cap: usize) -> usize {
        let threads = threads.max(1) as u64;
        let per_thread = self.0 / (threads * BYTES_PER_SCORE);
        (per_thread as usize).clamp(MIN_CHUNK_SCORES, cap)
    }
}

/// Never read a window smaller than this many bases; tiny windows re-read overlap out of
/// proportion to the sequence they cover.
const MIN_WINDOW_BASES: usize = 1 << 20;

impl MemoryBudget {
    /// How many bases of sequence to hold at once when reading through an index.
    ///
    /// The window is the sequence side of the footprint, so it gets a share of the budget
    /// separate from the score buffers, and is capped because past a point a larger
    /// window only costs memory without saving reads.
    pub fn window_bases(&self, cap: usize) -> usize {
        // A quarter of the budget: the scores buffered from a window dominate it, at
        // eight bytes per base against the one byte the base itself takes.
        let share = (self.0 / 4) as usize;
        share.clamp(MIN_WINDOW_BASES, cap)
    }
}

impl fmt::Display for MemoryBudget {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let b = self.0;
        if b >= GI && b.is_multiple_of(GI) {
            write!(f, "{}G", b / GI)
        } else if b >= MI && b.is_multiple_of(MI) {
            write!(f, "{}M", b / MI)
        } else if b >= KI && b.is_multiple_of(KI) {
            write!(f, "{}K", b / KI)
        } else {
            write!(f, "{b}")
        }
    }
}

/// The string was not a memory size this tool understands.
#[derive(Debug, PartialEq, Eq)]
pub struct ParseBudgetError(String);

impl fmt::Display for ParseBudgetError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}: expected a size like 8G, 512M, 64K or a plain byte count",
            self.0
        )
    }
}

impl std::error::Error for ParseBudgetError {}

impl FromStr for MemoryBudget {
    type Err = ParseBudgetError;

    /// Parse a size such as `8G`, `512M`, `64K` or `1048576`.
    ///
    /// Suffixes are binary: `8G` is 8 * 1024^3. `KB`/`KiB` spellings are accepted and
    /// mean the same thing, since offering both binary and decimal units under names
    /// that differ by one letter invites mistakes.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let err = || ParseBudgetError(s.to_string());
        let t = s.trim();
        if t.is_empty() {
            return Err(err());
        }
        // Strip an optional trailing "b"/"ib" so KB, KiB and K all parse alike.
        let t = t.strip_suffix(['b', 'B']).unwrap_or(t);
        let t = t.strip_suffix(['i', 'I']).unwrap_or(t);

        let (digits, multiplier) = match t.chars().last().ok_or_else(err)? {
            'k' | 'K' => (&t[..t.len() - 1], KI),
            'm' | 'M' => (&t[..t.len() - 1], MI),
            'g' | 'G' => (&t[..t.len() - 1], GI),
            't' | 'T' => (&t[..t.len() - 1], TI),
            _ => (t, 1),
        };

        let value: u64 = digits.trim().parse().map_err(|_| err())?;
        let bytes = value.checked_mul(multiplier).ok_or_else(err)?;
        if bytes == 0 {
            return Err(err());
        }
        Ok(MemoryBudget(bytes))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parses_suffixes() {
        let cases = [
            ("8G", 8 * GI),
            ("8g", 8 * GI),
            ("8GB", 8 * GI),
            ("8GiB", 8 * GI),
            ("512M", 512 * MI),
            ("64K", 64 * KI),
            ("2T", 2 * TI),
            ("1048576", 1048576),
            ("  4G  ", 4 * GI),
        ];
        for (text, expected) in cases {
            assert_eq!(
                text.parse::<MemoryBudget>().map(|b| b.bytes()),
                Ok(expected),
                "for {text:?}"
            );
        }
    }

    #[test]
    fn test_rejects_nonsense() {
        for text in ["", "   ", "G", "8X", "-1", "1.5G", "eight", "0", "0G"] {
            assert!(
                text.parse::<MemoryBudget>().is_err(),
                "{text:?} should not parse"
            );
        }
    }

    #[test]
    fn test_rejects_overflow_rather_than_wrapping() {
        assert!("99999999999T".parse::<MemoryBudget>().is_err());
    }

    #[test]
    fn test_display_round_trips() {
        for text in ["8G", "512M", "64K", "1023"] {
            let parsed: MemoryBudget = text.parse().unwrap();
            assert_eq!(parsed.to_string(), text);
        }
    }

    #[test]
    fn test_chunk_scores_scales_with_budget_and_threads() {
        let cap = 1 << 20;
        let eight_g: MemoryBudget = "8G".parse().unwrap();
        // 8 GiB over 10 threads is far more than the cap, so the cap applies.
        assert_eq!(eight_g.chunk_scores(10, cap), cap);

        // A small budget actually binds: 64 MiB / (8 threads * 8 bytes) = 1 Mi scores,
        // which is just under the cap.
        let small: MemoryBudget = "64M".parse().unwrap();
        assert_eq!(small.chunk_scores(8, cap), (64 * MI / (8 * 8)) as usize);

        // A tiny budget floors rather than collapsing to nothing.
        let tiny: MemoryBudget = "1K".parse().unwrap();
        assert_eq!(tiny.chunk_scores(16, cap), MIN_CHUNK_SCORES);
    }

    #[test]
    fn test_window_bases_scales_and_clamps() {
        let cap = 64 << 20;
        let big: MemoryBudget = "8G".parse().unwrap();
        assert_eq!(big.window_bases(cap), cap, "a large budget hits the cap");

        let mid: MemoryBudget = "64M".parse().unwrap();
        assert_eq!(mid.window_bases(cap), (64 * MI / 4) as usize);

        let tiny: MemoryBudget = "1K".parse().unwrap();
        assert_eq!(
            tiny.window_bases(cap),
            MIN_WINDOW_BASES,
            "floors rather than collapsing"
        );
    }

    #[test]
    fn test_chunk_scores_handles_zero_threads() {
        let b: MemoryBudget = "8G".parse().unwrap();
        assert_eq!(b.chunk_scores(0, 1 << 20), 1 << 20);
    }
}
