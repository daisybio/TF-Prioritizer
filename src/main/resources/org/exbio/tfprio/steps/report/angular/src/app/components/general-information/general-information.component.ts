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
}
