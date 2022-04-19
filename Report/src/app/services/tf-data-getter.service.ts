import {Injectable} from '@angular/core';
import {InputData, TranscriptionFactorGroup} from "../types/types";
//import data from '../../assets/data.json';

import data from '../../assets/input/data.json';


@Injectable({
  providedIn: 'root'
})
export class TfDataGetterService {
  inputData: InputData;
  tfMap: { [name: string]: TranscriptionFactorGroup } = {};

  constructor() {
    // @ts-ignore
    this.inputData = data;

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
