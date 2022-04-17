import {Injectable} from '@angular/core';
import {Data} from "../types/types";

@Injectable({
  providedIn: 'root'
})
export class MapToTableService {

  constructor() {
  }

  mapToTable(data: Data[]) {
    let matrix: number[][];
    let accessors: string[] = [];

    const dimensions = data[0].groups.length;
    if (dimensions < 1 || dimensions > 2) {
      console.error("Cannot process data with " + dimensions + " dimensions");
    }

    for (let entry of data) {
      if (entry.groups.length != dimensions) {
        console.error("Not all date entries have same dimension number");
      }

      accessors.push(...entry.groups);
    }
    accessors.sort();

    for (let entry of data) {

    }

    return accessors;
  }
}
