import {Component} from '@angular/core';
import {DataManagerService} from "../../../services/data-manager.service";

@Component({
  selector: 'app-parameters',
  templateUrl: './parameters.component.html',
  styleUrls: ['./parameters.component.scss']
})
export class ParametersComponent {
  configs = this.dataService.getConfigs();

  constructor(private dataService: DataManagerService) {
  }

  getEntries(value: any): any[] {
    return Object.entries(value);
  }

  isObject(value: any): boolean {
    return typeof value === 'object' && value !== null;
  }
}
