import {Component, Input, OnInit} from '@angular/core';
import {TranscriptionFactorGroup} from "../types/types";
import {TfDataGetterService} from "../services/tf-data-getter.service";

@Component({
  selector: 'tf-ranking',
  templateUrl: './tf-ranking.component.html',
  styleUrls: ['./tf-ranking.component.css']
})
export class TfRankingComponent implements OnInit {
  @Input()
  tfGroups: TranscriptionFactorGroup[];

  constructor(tfDataGetter: TfDataGetterService) {
    this.tfGroups = tfDataGetter.getData().transcriptionFactorGroups;
  }

  ngOnInit(): void {
  }
}
