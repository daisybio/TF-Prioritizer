import {Component, OnInit} from '@angular/core';
import {TranscriptionFactorGroup} from "../../types/types";
import {ActivatedRoute} from "@angular/router";
import {TfDataGetterService} from "../../services/tf-data-getter.service";

@Component({
  selector: 'app-regression',
  templateUrl: './regression.component.html',
  styleUrls: ['./regression.component.css']
})
export class RegressionComponent implements OnInit {
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
