import { Component } from '@angular/core';
import {DataManagerService} from "./services/data-manager.service";
import {tfGroup} from "../interfaces";

@Component({
  selector: 'app-root',
  templateUrl: './app.component.html',
  styleUrls: ['./app.component.sass']
})
export class AppComponent {
  groups: tfGroup[] = this.dataService.groups;

  constructor(private dataService: DataManagerService) {
  }
}
