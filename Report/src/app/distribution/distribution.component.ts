import {Component, OnInit} from '@angular/core';
import {TranscriptionFactorGroup} from "../types/types";
import {ActivatedRoute} from "@angular/router";
import {TfDataGetterService} from "../services/tf-data-getter.service";
import {InformationGetterService} from "../services/information-getter.service";

@Component({
  selector: 'app-distribution',
  templateUrl: './distribution.component.html',
  styleUrls: ['./distribution.component.css']
})

export class DistributionComponent implements OnInit {
  tfGroup: TranscriptionFactorGroup | undefined;
  routedName: string | null;
  information;

  dcgVisible: boolean = false;

  histoneModifications: string[] = [];

  constructor(private route: ActivatedRoute, private tfGetter: TfDataGetterService, private informationGetter: InformationGetterService) {
    this.routedName = this.route.snapshot.queryParamMap.get("tf");
    if (this.routedName) {
      this.tfGroup = tfGetter.getTfGroupByName(this.routedName);
    }
    this.information = informationGetter.getInformation()["distribution"];

    if (this.tfGroup) {
      Object.keys(this.tfGroup.distribution.plots).forEach(hm => this.histoneModifications.push(hm));
    }
  }

  getRank(hm: string) {
    return this.tfGroup?.distribution.ranks[hm]["rank"];
  }

  getSize(hm: string) {
    return this.tfGroup?.distribution.ranks[hm]["size"];
  }

  ngOnInit(): void {
  }
}
