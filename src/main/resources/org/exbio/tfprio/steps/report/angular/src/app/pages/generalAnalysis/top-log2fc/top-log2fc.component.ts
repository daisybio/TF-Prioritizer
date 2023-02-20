import {Component} from '@angular/core';
import {DataManagerService} from "../../../services/data-manager.service";

@Component({
  selector: 'app-top-log2fc',
  templateUrl: './top-log2fc.component.html',
  styleUrls: ['./top-log2fc.component.scss']
})
export class TopLog2fcComponent {
  topLog2fc;

  public constructor(dataService: DataManagerService) {
    this.topLog2fc = dataService.getTopLog2fc();
  }
}
