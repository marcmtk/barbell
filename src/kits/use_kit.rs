use crate::annotate::annotator::annotate_with_kit;
use crate::config::{AnnotateConfig, FilterConfig, KitConfig, TrimConfig};
use crate::filter::filter::filter;
use crate::inspect::inspect::inspect;
use crate::kits::kits::*;
use crate::trim::trim::{LabelSide, trim_matches};
use anyhow::anyhow;
use colored::*;
use needletail::parse_fastx_file;
use pa_types::Cost;
use sassy::Searcher;
use sassy::profiles::Iupac;
use std::path::{Path, PathBuf};

// Reads, edits and minimum fraction used to detect the SQK-MAB114.24 amplicon
const MAB_DETECT_READS: usize = 10_000;
const MAB_DETECT_MAX_EDITS: usize = 8;
const MAB_DETECT_MIN_FRAC: f64 = 0.95;

/// SQK-MAB114.24 is used for 16S and ITS amplicons, detect which one from the primer
/// that follows the barcode cassette. Other kits are returned unchanged
pub fn resolve_kit_name(kit: &str, fastq_files: &[PathBuf]) -> anyhow::Result<String> {
    if kit.replace('.', "-") != "SQK-MAB114-24" {
        return Ok(kit.to_string());
    }
    let (n_16s, n_its) = count_mab_amplicons(fastq_files)?;
    let total = (n_16s + n_its) as f64;
    println!("Amplicons detected in the first {MAB_DETECT_READS} reads: 16S={n_16s}, ITS={n_its}");
    let resolved = if total >= 100.0 && n_16s as f64 >= MAB_DETECT_MIN_FRAC * total {
        "SQK-MAB114-24-16S"
    } else if total >= 100.0 && n_its as f64 >= MAB_DETECT_MIN_FRAC * total {
        "SQK-MAB114-24-ITS"
    } else {
        println!(
            "{}",
            "Could not determine 16S or ITS, concatenated reads will not be split. \
             Use SQK-MAB114-24-16S or SQK-MAB114-24-ITS to set it"
                .yellow()
        );
        "SQK-MAB114-24"
    };
    println!("Using kit: {resolved}");
    Ok(resolved.to_string())
}

fn count_mab_amplicons(fastq_files: &[PathBuf]) -> anyhow::Result<(usize, usize)> {
    let mut searcher = Searcher::<Iupac>::new_fwd();
    let (mut n_reads, mut n_16s, mut n_its) = (0, 0, 0);
    for path in fastq_files {
        let mut reader = parse_fastx_file(path)?;
        while let Some(record) = reader.next() {
            if n_reads == MAB_DETECT_READS {
                return Ok((n_16s, n_its));
            }
            n_reads += 1;
            let record = record?;
            let seq = record.seq();
            let start = &seq[..seq.len().min(250)];
            let cost_16s = best_cost(&mut searcher, &[MAB_REAR_16S_F, MAB_REAR_16S_R], start);
            let cost_its = best_cost(&mut searcher, &[MAB_REAR_ITS1, MAB_REAR_ITS4], start);
            match (cost_16s, cost_its) {
                (Some(a), Some(b)) if a < b => n_16s += 1,
                (Some(a), Some(b)) if b < a => n_its += 1,
                (Some(_), None) => n_16s += 1,
                (None, Some(_)) => n_its += 1,
                _ => {}
            }
        }
    }
    Ok((n_16s, n_its))
}

fn best_cost(searcher: &mut Searcher<Iupac>, flanks: &[&str], text: &[u8]) -> Option<Cost> {
    flanks
        .iter()
        .flat_map(|flank| searcher.search(flank.as_bytes(), text, MAB_DETECT_MAX_EDITS))
        .map(|m| m.cost)
        .min()
}

pub fn demux_using_kit(fastq_files: &[PathBuf], config: &KitConfig) -> anyhow::Result<()> {
    let kit_name = resolve_kit_name(config.kit_name.as_str(), fastq_files)?;
    let kit_name = kit_name.as_str();
    let output_folder = config.output_folder.as_str();
    // Create output folder if not exists yet
    if !Path::new(output_folder).exists() {
        std::fs::create_dir_all(output_folder)?;
    }

    let kit_info = get_kit_info(kit_name);

    // Print some kit info
    println!("\n{}", "Kit info".purple().bold());
    println!("Kit name: {}", kit_info.name);
    println!(
        "Kit type: {}",
        if config.maximize { "Maximize" } else { "Safe" }
    );
    for tmpl in kit_info.templates {
        println!("Barcodes: {} - {}", tmpl.barcodes.from, tmpl.barcodes.to);
    }

    // If the default values are
    println!("\n{}", "Annotating reads...".purple().bold());
    let annotate_config = AnnotateConfig {
        max_flank_errors: config.max_flank_errors,
        alpha: config.alpha,
        n_threads: config.threads as u32,
        verbose: config.verbose,
        min_score: config.min_score,
        min_score_diff: config.min_score_diff,
        use_extended: config.use_extended,
    };
    annotate_with_kit(
        fastq_files,
        format!("{output_folder}/annotation.tsv").as_str(),
        kit_name,
        &annotate_config,
    )?;

    // // After annotating we show inspect
    println!("\n{}", "Top 10 most common patterns".purple().bold());
    let pattern_per_read_out = format!("{output_folder}/pattern_per_read.tsv");
    inspect(
        format!("{output_folder}/annotation.tsv").as_str(),
        10,
        Some(pattern_per_read_out),
        250,
    )
    .map_err(|e| anyhow!("{e}"))?;
    println!(
        "Want to see more patterns? Run: `barbell inspect -i {output_folder}/annotation.tsv -n 100`"
    );

    // Filter
    println!("\n{}", "Filtering reads...".purple().bold());

    let patterns = if config.maximize {
        (kit_info.maximize_patterns)()
    } else {
        (kit_info.safe_patterns)()
    };
    let filter_config = FilterConfig {
        verbose: config.verbose,
    };

    filter(
        format!("{output_folder}/annotation.tsv").as_str(),
        format!("{output_folder}/filtered.tsv").as_str(),
        None,
        patterns,
        &filter_config,
    )
    .map_err(|e| anyhow!("{e}"))?;

    // Trimming
    println!("\n{}", "Trimming reads...".purple().bold());
    let trim_config = TrimConfig {
        add_labels: true,
        add_orientation: false,
        add_flank: false,
        sort_labels: false,
        only_side: Some(LabelSide::Left),
        failed_trimmed_writer: config.failed_out.clone(),
        write_full_header: true,
        skip_trim: false,
        flip: false,
        verbose: config.verbose,
        gzip: config.gzip,
    };
    trim_matches(
        format!("{output_folder}/filtered.tsv").as_str(),
        fastq_files,
        output_folder,
        &trim_config,
    )?;

    println!("\n{}", "Done!".green().bold());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    // Barcode cassette (AB05) followed by a concrete primer, like the start of a MAB read
    fn mab_reads(primers: &[&str]) -> NamedTempFile {
        let mut tmp = NamedTempFile::with_suffix(".fastq").unwrap();
        for (i, primer) in primers.iter().enumerate() {
            let seq = format!(
                "ACGTACGTACGTACGTACGTTTTAACCAGACTCTTGGGCTCGCCACGTCCCCCATATCCGTGTCGCCCTT{primer}{}",
                "ACGTTGCA".repeat(25)
            );
            writeln!(tmp, "@read{i}\n{seq}\n+\n{}", "5".repeat(seq.len())).unwrap();
        }
        tmp.flush().unwrap();
        tmp
    }

    #[test]
    fn resolve_mab_amplicon() {
        let primer_16s = "AGAGTTTGATCATGGCTCAG"; // 16S_Ent_F
        let primer_its = "TCCGTAGGTGAACCTGCGG"; // ITS1
        let resolve = |primers: Vec<&str>| {
            let tmp = mab_reads(&primers);
            resolve_kit_name("SQK-MAB114.24", &[tmp.path().to_path_buf()]).unwrap()
        };
        assert_eq!(resolve(vec![primer_16s; 200]), "SQK-MAB114-24-16S");
        assert_eq!(resolve(vec![primer_its; 200]), "SQK-MAB114-24-ITS");
        // Mixed runs fall back to the combined kit
        let mixed = [vec![primer_16s; 100], vec![primer_its; 100]].concat();
        assert_eq!(resolve(mixed), "SQK-MAB114-24");
        // Too few reads to tell
        assert_eq!(resolve(vec![primer_16s; 20]), "SQK-MAB114-24");
        // Other kits are left alone
        assert_eq!(
            resolve_kit_name("SQK-NBD114-96", &[]).unwrap(),
            "SQK-NBD114-96"
        );
    }
}
