import {Injectable} from '@angular/core';
import data from '../../assets/data.json';
import {dataInterface, tfGroup} from "../../interfaces";

@Injectable({
  providedIn: 'root'
})
export class DataManagerService {
  // @ts-ignore
  private _inputData: dataInterface = data;

  private groupMap = this._inputData.groups.reduce((acc, group) => {
    acc.set(group.symbol, group);
    return acc;
  }, new Map<string, tfGroup>());

  constructor() {
  }

  get groups(): tfGroup[] {
    return this._inputData.groups;
  }

  get hms(): string[] {
    return Array.from(this.groups.flatMap(group => Object.keys(group.hmDcg))
      .reduce((acc, hm) => acc.add(hm), new Set<string>())).sort();
  }

  get regressionCoefficients() {
    return this._inputData.regressionCoefficients;
  }

  public getTfGroupByName(name: string): tfGroup | undefined {
    return this.groupMap.get(name);
  }

  public formatPlotData(data: { [p: string]: number }) {
    return Object.entries(data).map(entry => {
      return {"name": entry[0], "value": entry[1]}
    });
  }

  public getConfigs() {
    return this._inputData.configs;
  }
}
