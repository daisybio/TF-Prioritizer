import {Component, OnInit} from '@angular/core';
import {TranscriptionFactorGroup} from "../../types/types";
import {ActivatedRoute} from "@angular/router";
import {TfDataGetterService} from "../../services/tf-data-getter.service";
import {InformationGetterService} from "../../services/information-getter.service";
import {formatNumber} from "@angular/common";

@Component({
  selector: 'app-regression',
  templateUrl: './regression.component.html',
  styleUrls: ['./regression.component.css']
})
export class RegressionComponent implements OnInit {
  tfGroup!: TranscriptionFactorGroup;
  routedName: string | null;
  information;
  tableVisible: boolean = false;
  histoneModifications: string[] | undefined;
  groupPairings: string[] | undefined;

  constructor(private route: ActivatedRoute, private tfGetter: TfDataGetterService, private informationGetter: InformationGetterService) {
    this.routedName = this.route.snapshot.queryParamMap.get("tf");
    if (this.routedName) {
      this.tfGroup = tfGetter.getTfGroupByName(this.routedName);
      this.histoneModifications = Object.keys(this.tfGroup.regression.table)
      let pairingSet = new Set<string>();
      this.histoneModifications.forEach(hm => Object.keys(this.tfGroup.regression.table[hm]).forEach(pairing => pairingSet.add(pairing)))
      this.groupPairings = Array.from<string>(pairingSet.values()).sort();
    }
    this.information = informationGetter.getInformation()["regression"];
  }

  ngOnInit(): void {
  }

  getValue(hm: string, groupPairing: string) {
    if (!this.tfGroup.regression.table[hm][groupPairing]) {
      return undefined;
    }
    return formatNumber(this.tfGroup.regression.table[hm][groupPairing], "en-GB", ".2");
  }
}
