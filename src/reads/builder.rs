use std::io::Write;

use crate::writer::AlignmentWriter;

#[derive(Debug, Clone, Copy)]
pub enum Flag {
    #[allow(dead_code)]
    Unmapped = 0x4,
    ReverseComplement = 0x10,
    SecondaryAlignment = 0x100,
    SupplementaryAlignment = 0x800,
}

pub struct TagAndValue {
    tag: String,
    value: TagValue,
}

impl<T: Into<TagValue>> From<(String, T)> for TagAndValue {
    fn from((tag, value): (String, T)) -> Self {
        TagAndValue {
            tag,
            value: value.into(),
        }
    }
}

impl<T: Into<TagValue>> From<(&str, T)> for TagAndValue {
    fn from((tag, value): (&str, T)) -> Self {
        TagAndValue {
            tag: tag.to_string(),
            value: value.into(),
        }
    }
}

pub enum TagValue {
    String(String),
    Integer(i64),
    Float(f64),
}

impl From<i32> for TagValue {
    fn from(value: i32) -> Self {
        TagValue::Integer(value as i64)
    }
}

impl From<i64> for TagValue {
    fn from(value: i64) -> Self {
        TagValue::Integer(value)
    }
}

impl From<usize> for TagValue {
    fn from(value: usize) -> Self {
        TagValue::Integer(value as i64)
    }
}

impl From<f32> for TagValue {
    fn from(value: f32) -> Self {
        TagValue::Float(value as f64)
    }
}

impl From<String> for TagValue {
    fn from(value: String) -> Self {
        TagValue::String(value)
    }
}

impl From<&str> for TagValue {
    fn from(value: &str) -> Self {
        TagValue::String(value.to_string())
    }
}

impl std::fmt::Display for TagAndValue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self.value {
            TagValue::String(s) => write!(f, "{}:Z:{}", self.tag, s),
            TagValue::Integer(i) => write!(f, "{}:i:{}", self.tag, i),
            TagValue::Float(fl) => write!(f, "{}:f:{}", self.tag, fl),
        }
    }
}

pub struct SegmentBuilder<'a> {
    qname: &'a str,
    flag: u16,
    rname: Option<&'a str>,
    pos: Option<usize>,
    mapq: Option<u8>,
    cigar: Option<&'a str>,
    seq: Option<&'a [u8]>,
    qual: Option<&'a [u8]>,
    primary: Option<(&'a str, usize)>,
    tags: Vec<TagAndValue>,
}

impl<'a> SegmentBuilder<'a> {
    pub fn new(qname: &'a str) -> Self {
        Self {
            qname,
            flag: 0,
            rname: None,
            pos: None,
            mapq: None,
            cigar: None,
            seq: None,
            qual: None,
            primary: None,
            tags: Vec::new(),
        }
    }

    pub fn set_flag(&mut self, flag: Flag) {
        self.flag |= flag as u16;
    }

    pub fn with_flags(mut self, flags: &[Flag]) -> Self {
        for flag in flags {
            self.set_flag(*flag);
        }
        self
    }

    pub fn with_reference(mut self, rname: &'a str, pos: usize) -> Self {
        self.rname = Some(rname);
        self.pos = Some(pos);
        self
    }

    pub fn with_mapping_quality(mut self, mapq: u8) -> Self {
        self.mapq = Some(mapq);
        self
    }

    pub fn with_cigar(mut self, cigar: &'a str) -> Self {
        self.cigar = Some(cigar);
        self
    }

    pub fn with_primary(mut self, primary: Option<(&'a str, usize)>) -> Self {
        self.primary = primary;
        self
    }

    pub fn with_sequence_and_quality(mut self, seq: &'a [u8], qual: &'a [u8]) -> Self {
        self.seq = Some(seq);
        self.qual = Some(qual);
        self
    }

    pub fn with_tag_and_value<T: Into<TagValue>>(mut self, tag: &str, value: T) -> Self {
        self.tags.push(TagAndValue::from((tag, value)));
        self
    }

    pub fn write<W: Write>(self, writer: &AlignmentWriter<W>) -> std::io::Result<()> {
        let rname = self.rname.unwrap_or("*");
        let pos = self.pos.unwrap_or(0);
        let mapq = self.mapq.unwrap_or(255);
        let cigar_str = self.cigar.unwrap_or("*");
        let (rnext, pnext) = if let Some((primary_rname, primary_pos)) = self.primary {
            (primary_rname, primary_pos)
        } else {
            ("*", 0)
        };
        let seq_str = if let Some(seq) = self.seq {
            String::from_utf8_lossy(seq).to_string()
        } else {
            "*".to_string()
        };
        let qual_str = if let Some(qual) = self.qual {
            String::from_utf8_lossy(qual).to_string()
        } else {
            "*".to_string()
        };
        let tags_str = self
            .tags
            .iter()
            .map(|tag| tag.to_string())
            .collect::<Vec<_>>()
            .join("\t");
        writer.write_alignment(
            self.qname,
            self.flag,
            rname,
            pos,
            mapq,
            &cigar_str,
            rnext,
            pnext,
            0, // tlen is not calculated here
            &seq_str,
            &qual_str,
            &tags_str,
        )?;
        Ok(())
    }
}
