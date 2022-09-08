import {Injectable} from '@angular/core';
import {InputData, TranscriptionFactorGroup} from "../types/types";

import {tfData} from "../../assets/input/data";


@Injectable({
  providedIn: 'root'
})
export class TfDataGetterService {
  tfMap: { [name: string]: TranscriptionFactorGroup } = {};
  private readonly inputData: InputData;

  constructor() {
    // @ts-ignore
    this.inputData = tfData;

    for (let transcriptionFactorGroup of this.inputData.transcriptionFactorGroups) {
      this.tfMap[transcriptionFactorGroup.name] = transcriptionFactorGroup;
    }
  }

  getData() {
    return this.inputData;
  }

  getTfGroupByName(name: string) {
    return this.tfMap[name];
  }

  getRanked(hms: string[]) {
    let name_score: {
      [name: string]: number
    } = {}

    for (let hm of hms) {
      let sorted = this.getRankedForHm(hm);
      console.log(hm);
      console.log(sorted);

      for (let i = 0; i < sorted.length; i++) {
        let tf = sorted[i];
        if (!Object.keys(name_score).includes(tf)) name_score[tf] = 0;

        name_score[tf] = name_score[tf] + (sorted.length - i) / sorted.length;
      }
    }

    let sortable: (string | number)[][] = [];
    for (let tf in name_score) {
      sortable.push([tf, name_score[tf]]);
    }

    sortable.sort((a, b) => {
      // @ts-ignore
      return a[1] - b[1];
    });
    // @ts-ignore
    let sorted: string[] = sortable.map(x => x[0]);
    return sorted.reverse();
  }

  getRankedForHm(hm: string) {
    let allTfs = Object.keys(this.tfMap);
    let withData: string[] = [];

    for (let tf of allTfs) {
      if (Object.keys(this.getTfGroupByName(tf).stats).includes(hm)) {
        withData.push(tf);
      }
    }

    return withData.sort((a: string, b: string) => {
      let aStats = this.getTfGroupByName(a).stats[hm];
      let bStats = this.getTfGroupByName(b).stats[hm];

      let medianDifference = aStats.median - bStats.median;
      if (medianDifference != 0) return medianDifference;

      let quantile95Difference = aStats["95_quantile"] - bStats["95_quantile"];
      if (quantile95Difference != 0) return quantile95Difference;

      return aStats["99_quantile"] - bStats["99_quantile"];
    }).reverse();
  }
}
