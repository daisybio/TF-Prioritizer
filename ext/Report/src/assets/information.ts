export const loremIpsum = "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Duis a lorem neque. Integer condimentum purus a maximus vehicula. Aenean interdum ligula vitae interdum ullamcorper. Sed vehicula sem id leo condimentum elementum. Ut varius quis turpis vitae tristique. Donec viverra venenatis condimentum. Aliquam sit amet metus nisl.";

export const INFORMATION = {
  "basicData": {
    "Log2fc": {
      "content": "The log-2-fold-changes of the given transcription factor for each group pairing calculated by DESeq2.",
      "source": ""
    },
    "NormEx": {
      "content": "The normalized transcription factor expression is calculated by DeSeq2.",
      "source": ""
    },
    "tpm": {
      "content": "Transcripts per million. Describes the transcription rate of the given transcription factor.",
      "source": ""
    }
  },
  "validation": {
    "logos": {
      "general": {
        "content": "The sequence logos of the transcription factor binding sequence.",
        "source": ""
      },
      "biophysical": {
        "content": "The sequence logos, which are calculated based on the TEPIC position weight matrix files",
        "source": ""
      },
      "tfSequence": {
        "content": "The binding sequence logo obtained from JASPAR.",
        "source": "Castro-Mondragon JA, Riudavets-Puig R, Rauluseviciute I, Berhanu Lemma R, Turchi L, Blanc-Mathieu R, Lucas J, Boddie P, Khan A, Manosalva Pérez N, Fornes O, Leung TY, Aguirre A, Hammal F, Schmelter D, Baranasic D, Ballester B, Sandelin A, Lenhard B, Vandepoele K, Wasserman WW, Parcy F, and Mathelier A JASPAR 2022: the 9th release of the open-access database of transcription factor binding profiles Nucleic Acids Res. 2022 Jan 7;50(D1):D165-D173.; doi: 10.1093/nar/gkab1113"
      }
    },
    "heatmap": {
      "content": "These heatmaps display the batch-corrected differential expression of the [TF]’s top target genes. Their goal is to visualize clusterings of the compared time points.",
      "source": ""
    },
    "igv": {
      "content": "TF-Prioritizer allows users to manually investigate the ChIP-seq signal in the identified CREs of differentially expressed genes. TF-prioritizer generates a compendium of screenshots of the top 30 upregulated or downregulated loci (sorted by their total log fold change) between two sample groups.",
      "source": ""
    }
  },
  "distribution": {
    "plots": {
      "content": "To determine which TFs have a significant contribution to a condition-specific change between two sample groups, we want to consider multiple lines of evidence in an aggregated score. We introduce Transcription Factor Target Gene scores (TF-TG scores, Figure 2) which combine (i) the absolute log2 fold change of differentially expressed genes since genes showing large expression differences are more likely affected through TF regulation than genes showing only minor expression differences; (ii), the TF-Gene scores from TEPIC indicating which TFs likely influence a gene, and (iii) to further quantify this link we also consider the total coefficients of a logistic regression model computed with DYNAMITE.",
      "source": ""
    },
    "ranks": {
      "content": "In order to achieve one list of prioritized TFs globally ranked over several histone modifications, a discounted cumulative gain (DCG) approach is being used. Here, the rank of the TF within each distinct HM is shown.",
      "source": ""
    }
  },
  "regression": {
    "coefficients": {
      "content": loremIpsum,
      "source": ""
    },
    "heatmaps": {
      "content": loremIpsum,
      "source": ""
    },
    "table": {
      "content": loremIpsum,
      "source": ""
    }
  }
}
