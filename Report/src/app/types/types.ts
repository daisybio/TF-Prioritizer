export interface InputData {
  "transcriptionFactorGroups": TranscriptionFactorGroup[]
}

export interface TranscriptionFactorGroup {
  "transcriptionFactors": TranscriptionFactor[],
  "name": string,
  "targetGenes": TargetGene[],
  "validation": {
    "heatmap": {},
    "igv": {}
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
