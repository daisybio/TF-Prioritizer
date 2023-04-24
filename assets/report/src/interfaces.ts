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
