import {Component, OnInit} from '@angular/core';
import {ConfusionMatrix} from "../../types/types";
import {TfDataGetterService} from "../../services/tf-data-getter.service";

@Component({
  selector: 'app-statistical-evaluation',
  templateUrl: './statistical-evaluation.component.html',
  styleUrls: ['./statistical-evaluation.component.css']
})
export class StatisticalEvaluationComponent implements OnInit {
  statisticalEvaluation: {
    [tf: string]: {
      "chip"?: ConfusionMatrix,
      "experimental"?: ConfusionMatrix,
      "combined"?: ConfusionMatrix
    }
  }

  overviewVisible = true;

  private visibilities: {
    [tfName: string]: boolean
  } = {}

  constructor(dataGetter: TfDataGetterService) {
    this.statisticalEvaluation = dataGetter.getData().statisticalEvaluation;
    this.getTfs().forEach(tfName => this.visibilities[tfName] = false);
  }

  ngOnInit(): void {
  }

  getTfs() {
    let tfs = Object.keys(this.statisticalEvaluation).sort();
    tfs = tfs.filter(function (value, index, arr) {
      return value != "metricsPlot";
    })
    return tfs;
  }

  isVisible(tfName: string) {
    return this.visibilities[tfName];
  }

  toggleOverview() {
    this.overviewVisible = !this.overviewVisible;
  }

  toggleVisibility(tfName: string) {
    this.visibilities[tfName] = !this.visibilities[tfName];
  }

  getTitle(raw: string) {
    switch (raw) {
      case "chip":
        return "ChIP vs predicted";
      case "experimental":
        return "Experimental vs predicted";
      case "combined":
        return "ChIP and experimental vs predicted";
    }
    return "Unknown title"
  }
}
