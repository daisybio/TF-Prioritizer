import {Injectable} from '@angular/core'
import {InputData} from "../types/types";
import data from '../../assets/data.json';


@Injectable({
  providedIn: 'root'
})
export class TfDataGetterService {
  inputData: InputData;

  constructor() {
    // @ts-ignore
    this.inputData = data;
  }

  getData() {
    return this.inputData;
  }
}
