import {Component, Input, OnInit} from '@angular/core';
import {TfDataGetterService} from "../../services/tf-data-getter.service";

@Component({
  selector: 'app-parameters',
  templateUrl: './parameters.component.html',
  styleUrls: ['./parameters.component.css']
})
export class ParametersComponent implements OnInit {
  @Input()
  configs: {
    [module: string]: {
      [parameter: string]: boolean | number | [] | string
    }
  }

  modules: string[] = []

  visibilities: {
    [module: string]: boolean
  } = {}

  constructor(tfGetter: TfDataGetterService) {
    this.configs = tfGetter.getData().configs;

    for (let module of Object.keys(this.configs)) {
      this.modules.push(module);

      this.visibilities[module] = false;
    }

    this.modules.sort();
  }

  getParameters(module: string) {
    return Object.keys(this.configs[module]);
  }

  getValue(module: string, parameter: string) {
    return this.configs[module][parameter];
  }

  ngOnInit(): void {
  }
}
