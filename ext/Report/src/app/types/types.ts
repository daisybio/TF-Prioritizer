export interface InputData {
  "transcriptionFactorGroups": TranscriptionFactorGroup[],
  "configs": {
    [module: string]: {
      [parameter: string]: boolean | number | [] | string
    }
  },
  "importantLoci": {},
  "topLog2fc": {},
  "backgroundStats": Stats,
  "coOccurrenceAnalysis": Data[],
  "existingValues": {
    "groupPairing": string[],
    "hm": string[],
    "group": string[]
  },
  "statisticalEvaluation": {
    [tf: string]: {
      "chip"?: ConfusionMatrix,
      "experimental"?: ConfusionMatrix,
      "combined"?: ConfusionMatrix
    }
  }
}

export interface TranscriptionFactorGroup {
  "transcriptionFactors": TranscriptionFactor[],
  "name": string,
  "targetGenes": TargetGene[],
  "stats": Stats,
  "validation": {
    "heatmap": {},
    "igv": {},
    "logos": {
      "biophysical": {},
      "tfSequence": {}
    }
  },
  "distribution": {
    "plots": {},
    "ranks": {
      [hm: string]: {
        [entryType: string]: number
      }
    }
  },
  "regression": {
    "coefficients": {},
    "heatmaps": {},
    "table": {
      [hm: string]: {
        [groupPairing: string]: number
      }
    }
  }
}

export interface TranscriptionFactor {
  "name": string,
  "tpm": Data[],
  "normalizedExpression": Data[],
  "geneID": string,
  "log2fc": Data[]
}

export interface TargetGene {
  "symbol": string,
  "geneID": string
}

export interface Data {
  "groups": string[],
  "value": number
}

export interface Stats {
  [hm: string]: {
    "sum_all_values": number,
    "number_target_genes": number,
    "median": number,
    "99_quantile": number,
    "mean": number,
    "95_quantile": number
  }
}

export interface ConfusionMatrix {
  "fn": number,
  "tn": number,
  "fp": number,
  "tp": number
}
