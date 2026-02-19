use crate::{
    config,
    reads::seeds::{self, SeedCluster},
    utils::debug::{DebugFile, DebugOutput, DebugTsvWriter, TsvRow},
};

pub mod rmq_dp;
pub mod agglomerative;
pub mod kruskal;

// ── Debug file statics ──────────────────────────────────────────────────────

/// Debug TSV file with seeds grouped into clusters (before chaining).
static CLUSTERS_TSV: DebugFile<ClustersTsvDebug> = DebugFile::new();

// ── Concrete debug types ─────────────────────────────────────────────────────

type ClustersTsvRow<'a> = (&'a str, usize, usize, usize, usize, &'a str, usize, usize, &'a str, usize);

struct ClustersTsvDebug(DebugTsvWriter);

impl ClustersTsvDebug {
    const HEADERS: &[&str] = &[
        "read_name", "cluster_id", "read_start", "read_end", "read_len",
        "chrom", "ref_start", "ref_end", "strand", "match_len",
    ];
    const _CHECK: () = assert!(Self::HEADERS.len() == <ClustersTsvRow<'static> as TsvRow>::NUM_FIELDS);
}

impl DebugOutput for ClustersTsvDebug {
    type Item<'a> = ClustersTsvRow<'a>;
    fn create() -> Option<Self> {
        let _ = Self::_CHECK;
        let path = &config::get().seeding.debug_clusters_tsv;
        if path.is_empty() { return None; }
        let header = Self::HEADERS.join("\t");
        DebugTsvWriter::open(path, Some(&header)).ok().map(Self)
    }
    fn append(&self, item: &ClustersTsvRow<'_>) { self.0.append_row(item); }
    fn finish(&self) { self.0.finish(); }
}

pub fn write_clusters_debug(
    clusters: &[SeedCluster],
    read_name: &str,
    chrom_name: &str,
    strand_seq: &[u8],
    strand_qual: &[u8],
    read_len: usize,
    is_reverse: bool,
) {
    if false && seeds::CHAINS_SAM.is_enabled() {
        for (cluster_id, cluster) in clusters.iter().enumerate() {
            // Write debug chain SAM with SA tags linking seeds
            cluster.write_chain_sam(
                read_name,
                cluster_id, // cluster index as ID
                chrom_name,
                strand_seq,
                strand_qual,
            );
        }
    }

    // Write debug clusters TSV (seeds with cluster index)
    if CLUSTERS_TSV.is_enabled() {
        for (cluster_id, cluster) in clusters.iter().enumerate() {
            let strand = if is_reverse { "-" } else { "+" };
            for hit in cluster.chain.iter() {
                // Convert strand coordinates to forward coordinates
                let (fwd_start, fwd_end) = if is_reverse {
                    (read_len - hit.read_end(), read_len - hit.read_pos)
                } else {
                    (hit.read_pos, hit.read_end())
                };
                CLUSTERS_TSV.append(&(
                    read_name, cluster_id, fwd_start, fwd_end, strand_seq.len(),
                    chrom_name, hit.ref_pos, hit.ref_end(), strand, hit.match_len,
                ));
            }
        }
    }
}
