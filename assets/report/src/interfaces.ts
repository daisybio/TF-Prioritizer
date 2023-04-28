export interface ranking {
  [assay: string]: {
    [pairing: string]: string[]
  }
}

export interface dcg {
  [assay: string]: {
    [pairing: string]: number
  }
}

export interface tfData {
  heatmaps: {
    [pairing: string]: {
      [hm: string]: string
    }
  },
  model: string,
  jaspar_logos: {
    [matrix: string]: string
  },
  biophysical_logo: string,
  name: string,
  ensg: string
}

export interface expressionData {
  "TPM": {
    [condition: string]: number
  }
  "Counts": {
    [condition: string]: number
  }
  "DESeq": {
    [condition_pairing: string]: number
  }
}
