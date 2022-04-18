import {Component, Input, OnInit} from '@angular/core';
import {Data} from "../../../types/types";
import {formatNumber} from "@angular/common";

@Component({
  selector: 'data-content',
  templateUrl: './data-content.component.html',
  styleUrls: ['./data-content.component.css']
})
export class DataContentComponent implements OnInit {
  @Input()
    // @ts-ignore
  data: Data[];

  accessors: string[] = [];
  dimensions: number = 0;
  activeAccessors: string[] = [];

  containsAll = (first: string[], second: string[]) => second.every(entry => first.includes(entry));

  constructor() {
  }

  ngOnInit(): void {
    this.dimensions = this.data[0].groups.length;

    for (let entry of this.data) {
      if (entry.groups.length != this.dimensions) {
        console.error("Dimensions are not uniform")
      }
      this.accessors.push(...entry.groups);
    }
  }

  highLightAccessors(accessors: string[]) {
    if (accessors.length != this.dimensions && accessors.length != 0) {
      console.error("Accessors do not match dimensions");
    }

    this.activeAccessors = accessors;
  }

  getValue(accessors: string[]) {
    if (accessors.length != this.dimensions) {
      console.error("Accessors do not match dimensions");
    }

    let defaultValue = "-";

    if (accessors.length == 2 && accessors[0] == accessors[1]) {
      return defaultValue;
    }

    for (let entry of this.data) {
      if (this.containsAll(entry.groups, accessors)) {
        let value = entry.value;
        if (!value) {
          return defaultValue;
        }
        return formatNumber(entry.value, "en-GB", ".2");
      }
    }
    return defaultValue;
  }
}
