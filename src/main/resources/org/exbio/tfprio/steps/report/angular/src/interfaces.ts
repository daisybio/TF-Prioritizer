export interface tf {
  "symbol": string;
  "groupMeanExpression": {
    [group: string]: number;
  }

  ensgs: string[];
  "groupTpm": {
    [group: string]: number;
  }
  "pairingLog2fc": {
    [pairing: string]: number;
  }
}

export interface tfGroup {
  "symbol": string;
  "hmDcg": {
    [hm: string]: number;
  }
  "transcriptionFactors": tf[];
  "biophysicalLogo"?: string;
  "tfSequence": {};
  "regressionCoefficients": {
    [hm: string]: {
      [pairing: string]: number;
    }
  }
  "heatmaps": {};
  "igv": {};
  distributionPlots: {};
  confusionMatrix?: {
    truePositives: number;
    falsePositives: number;
    trueNegatives: number;
    falseNegatives: number;
  }
}

export interface dataInterface {
  "groups": tfGroup[];
  coOccurrence: {
    "name": string;
    "series": {
      "name": string;
      "value": number;
    }[];
  }[];
  "configs": {
    [mod: string]: {
      [config: string]: string | number | boolean;
    }
  }

  "importantLoci": {
    [group: string]: {
      [pattern: string]: {
        [locus: string]: string
      };
    }
  }

  "topLog2fc": {
    [pairing: string]: {
      "downregulated": {
        [gene: string]: string;
      };
      "upregulated": {
        [gene: string]: string;
      };
    }
  }
}

