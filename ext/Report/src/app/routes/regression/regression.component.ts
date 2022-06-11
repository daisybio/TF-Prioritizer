import {Component, OnInit} from '@angular/core';
import {TranscriptionFactorGroup} from "../../types/types";
import {ActivatedRoute} from "@angular/router";
import {TfDataGetterService} from "../../services/tf-data-getter.service";
import {InformationGetterService} from "../../services/information-getter.service";

@Component({
  selector: 'app-regression',
  templateUrl: './regression.component.html',
  styleUrls: ['./regression.component.css']
})
export class RegressionComponent implements OnInit {
  tfGroup: TranscriptionFactorGroup | undefined;
  routedName: string | null;
  information

  constructor(private route: ActivatedRoute, private tfGetter: TfDataGetterService, private informationGetter: InformationGetterService) {
    this.routedName = this.route.snapshot.queryParamMap.get("tf");
    if (this.routedName) {
      this.tfGroup = tfGetter.getTfGroupByName(this.routedName);
    }
    this.information = informationGetter.getInformation()["regression"];
  }

  ngOnInit(): void {
  }

}
