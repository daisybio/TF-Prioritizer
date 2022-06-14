import {Component, OnInit} from '@angular/core';
import {ActivatedRoute} from "@angular/router";
import {TranscriptionFactorGroup} from "../../types/types";
import {TfDataGetterService} from "../../services/tf-data-getter.service";
import {InformationGetterService} from "../../services/information-getter.service";

@Component({
  selector: 'app-validation',
  templateUrl: './validation.component.html',
  styleUrls: ['./validation.component.css']
})
export class ValidationComponent implements OnInit {
  tfGroup: TranscriptionFactorGroup | undefined;
  routedName: string | null;

  heatmapsVisible: boolean = false;
  igvVisible: boolean = false;
  logosVisible: boolean = false;
  logosBiophysicalVisible: boolean = false;
  logosTfSequenceVisible: boolean = false;
  information;

  constructor(private route: ActivatedRoute, private tfGetter: TfDataGetterService, private informationGetter: InformationGetterService) {
    this.routedName = this.route.snapshot.queryParamMap.get("tf");
    this.information = informationGetter.getInformation()['validation'];
    if (this.routedName) {
      this.tfGroup = tfGetter.getTfGroupByName(this.routedName);
    }
  }

  ngOnInit(): void {
  }
}
