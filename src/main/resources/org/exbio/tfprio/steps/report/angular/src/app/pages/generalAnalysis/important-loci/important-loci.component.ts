import {Component} from '@angular/core';
import {DataManagerService} from "../../../services/data-manager.service";

@Component({
  selector: 'app-important-loci',
  templateUrl: './important-loci.component.html',
  styleUrls: ['./important-loci.component.scss']
})
export class ImportantLociComponent {
  importantLoci;

  public constructor(dataService: DataManagerService) {
    this.importantLoci = dataService.getImportantLoci();
  }
}
