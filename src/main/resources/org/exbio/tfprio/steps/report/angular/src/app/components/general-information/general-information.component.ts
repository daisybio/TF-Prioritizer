import {ChangeDetectionStrategy, Component, Input} from '@angular/core';
import {tfGroup} from "../../../interfaces";
import {Subject} from "rxjs";
import {DataManagerService} from "../../services/data-manager.service";

@Component({
  selector: 'app-general-information',
  templateUrl: './general-information.component.html',
  styleUrls: ['./general-information.component.scss'],
  changeDetection: ChangeDetectionStrategy.OnPush
})
export class GeneralInformationComponent {
  @Input() tfGroup!: tfGroup;
  @Input() activeGroup?: Subject<tfGroup>;

  dataService!: DataManagerService;

  constructor(dataService: DataManagerService) {
    this.dataService = dataService;
  }

  isEmtpy(content: {}): boolean {
    // Check if any value is not 0
    return Object.values(content).every(value => value === 0);
  }

  getHeight(content: {}): {} {
    let size = Object.keys(content).length;
    let height = size * 50 + 30;
    return {
      'height': height + 'px'
    }
  }
}
