import {Component, OnInit} from '@angular/core';
import {TfDataGetterService} from "../../services/tf-data-getter.service";
import {Data} from "../../types/types";

@Component({
  selector: 'app-co-occurrence-analysis',
  templateUrl: './co-occurrence-analysis.component.html',
  styleUrls: ['./co-occurrence-analysis.component.css']
})
export class CoOccurrenceAnalysisComponent implements OnInit {
  data: Data[];

  constructor(private dataGetter: TfDataGetterService) {
    this.data = dataGetter.getData().coOccurrenceAnalysis;
    console.log(this.data.length);
  }

  ngOnInit(): void {
  }

}
