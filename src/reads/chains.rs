use crate::{
    reads::seeds::SeedCluster,
    utils::debug::{self, DebugFile},
};

pub mod rmq_dp;
pub mod agglomerative;
pub mod kruskal;

pub fn write_clusters_debug(
    clusters: &[SeedCluster],
    read_name: &str,
    chrom_name: &str,
    strand_seq: &[u8],
    strand_qual: &[u8],
    read_len: usize,
    is_reverse: bool,
) {
    if false && debug::is_enabled(DebugFile::ChainsSam) {
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
    if debug::is_enabled(DebugFile::ClustersTsv) {
        for (cluster_id, cluster) in clusters.iter().enumerate() {
            let strand = if is_reverse { "-" } else { "+" };
            for hit in cluster.chain.iter() {
                // Convert strand coordinates to forward coordinates
                let (fwd_start, fwd_end) = if is_reverse {
                    (read_len - hit.read_end(), read_len - hit.read_pos)
                } else {
                    (hit.read_pos, hit.read_end())
                };
                debug::write(
                    DebugFile::ClustersTsv,
                    &format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        read_name,
                        cluster_id,
                        fwd_start,
                        fwd_end,
                        strand_seq.len(),
                        chrom_name,
                        hit.ref_pos,
                        hit.ref_end(),
                        strand,
                        hit.match_len,
                    ),
                );
            }
        }
    }
}

