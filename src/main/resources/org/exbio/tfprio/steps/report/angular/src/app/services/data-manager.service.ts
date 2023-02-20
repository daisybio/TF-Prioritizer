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
    return this._inputData.groups.flatMap(group => Object.entries(group.regressionCoefficients));
  }

  get coOccurrence() {
    return this._inputData.coOccurrence;
  }

  public getTfGroupByName(name: string): tfGroup | undefined {
    return this.groupMap.get(name);
  }

  public formatPlotData(data: { [p: string]: number }) {
    return Object.entries(data).map(entry => {
      return {"name": entry[0], "value": entry[1]}
    });
  }

  public format2DPlotData(data: { [p: string]: { [p: string]: number } }) {
    return Object.entries(data).map(entry => {
      return {
        "name": entry[0],
        "series": this.formatPlotData(entry[1])
      }
    });
  }

  public getAllRegressionCoefficients() {
    let data: { [hm: string]: { [pairing: string]: { 'name': string, 'value': number }[] } } = {};

    this.groups.forEach(group => {
      let name = group.symbol;

      let regressionCoefficients = group.regressionCoefficients;

      Object.keys(regressionCoefficients).forEach(hm => {
        Object.keys(regressionCoefficients[hm]).forEach(pairing => {
          if (!data[hm]) {
            data[hm] = {};
          }

          if (!data[hm][pairing]) {
            data[hm][pairing] = [];
          }

          data[hm][pairing].push({
            "name": name,
            "value": regressionCoefficients[hm][pairing]
          });
        });
      });
    });

    return data;
  }

  public getConfigs() {
    return this._inputData.configs;
  }

  public getImportantLoci() {
    return this._inputData.importantLoci;
  }

  public getTopLog2fc() {
    return this._inputData.topLog2fc;
  }
}
