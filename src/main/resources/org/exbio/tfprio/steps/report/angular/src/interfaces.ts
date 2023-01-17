export interface tf {
  "symbol": string;
  "groupMeanExpression": {
    [group: string]: number;
  }

  ensgs: string[];
  "groupTpm": {
    [group: string]: number;
  }
}

export interface tfGroup {
  "symbol": string;
  "transcriptionFactors": tf[];
}

export interface dataInterface {
  "groups": tfGroup[];
}
