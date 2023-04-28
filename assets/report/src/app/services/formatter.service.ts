import {Injectable} from '@angular/core';

@Injectable({
  providedIn: 'root'
})
export class FormatterService {
  constructor() {
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
}
