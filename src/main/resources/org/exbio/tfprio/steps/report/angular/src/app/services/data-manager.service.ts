import { Injectable } from '@angular/core';
import data from '../../assets/data.json';
import {dataInterface, tfGroup} from "../../interfaces";
@Injectable({
  providedIn: 'root'
})
export class DataManagerService {
  private _inputData: dataInterface = data;

  get groups(): tfGroup[] {
    return this._inputData.groups;
  }

  constructor() { }
}
