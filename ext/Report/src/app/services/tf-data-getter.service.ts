import {Injectable} from '@angular/core';
import {InputData, TranscriptionFactorGroup} from "../types/types";

import {tfData} from "../../assets/input/data";


@Injectable({
  providedIn: 'root'
})
export class TfDataGetterService {
  private readonly inputData: InputData;
  tfMap: { [name: string]: TranscriptionFactorGroup } = {};

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
}
