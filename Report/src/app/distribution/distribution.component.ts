import {Component, OnInit} from '@angular/core';
import {TranscriptionFactorGroup} from "../types/types";
import {ActivatedRoute} from "@angular/router";
import {TfDataGetterService} from "../services/tf-data-getter.service";

@Component({
  selector: 'app-distribution',
  templateUrl: './distribution.component.html',
  styleUrls: ['./distribution.component.css']
})
export class DistributionComponent implements OnInit {
  tfGroup: TranscriptionFactorGroup | undefined;
  routedName: string | null;

  constructor(private route: ActivatedRoute, private tfGetter: TfDataGetterService) {
    this.routedName = this.route.snapshot.queryParamMap.get("tf");
    if (this.routedName) {
      this.tfGroup = tfGetter.getTfGroupByName(this.routedName);
    }
  }

  ngOnInit(): void {
  }

}
