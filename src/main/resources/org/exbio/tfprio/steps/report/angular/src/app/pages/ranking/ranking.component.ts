import { Component, OnInit } from '@angular/core';
import {tfGroup} from "../../../interfaces";
import {DataManagerService} from "../../services/data-manager.service";

@Component({
  selector: 'app-ranking',
  templateUrl: './ranking.component.html',
  styleUrls: ['./ranking.component.sass']
})
export class RankingComponent {
  groups: tfGroup[] = this.dataService.groups;

  constructor(private dataService: DataManagerService) {
  }

}
