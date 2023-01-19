import {Component} from '@angular/core';
import {ActivatedRoute} from "@angular/router";
import {DataManagerService} from "../../../services/data-manager.service";
import {tfGroup} from "../../../../interfaces";

@Component({
  selector: 'app-tf-details',
  templateUrl: './tf-details.component.html',
  styleUrls: ['./tf-details.component.scss']
})
export class TfDetailsComponent {
  tfGroup?: tfGroup;

  constructor(private route: ActivatedRoute, private dataManager: DataManagerService) {
    let routedName = this.route.snapshot.queryParamMap.get("tf");

    if (routedName) {
      this.tfGroup = dataManager.getTfGroupByName(routedName);
    }
  }
}
