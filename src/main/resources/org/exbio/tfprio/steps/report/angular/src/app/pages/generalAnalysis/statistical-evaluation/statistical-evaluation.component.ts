import {Component} from '@angular/core';
import {DataManagerService} from "../../../services/data-manager.service";
import {Subject} from "rxjs";

@Component({
  selector: 'app-statistical-evaluation',
  templateUrl: './statistical-evaluation.component.html',
  styleUrls: ['./statistical-evaluation.component.scss']
})
export class StatisticalEvaluationComponent {
  data: {
    [tf: string]: {
      truePositives: number;
      falsePositives: number;
      trueNegatives: number;
      falseNegatives: number;
    }
  }

  metrics: {
    [tf: string]: {
      Accuracy: number;
      Precision: number;
      Sensitivity: number;
      Specificity: number;
      F1_Score: number;
    }
  }

  $style: Subject<{}> = new Subject<{}>();

  plotData: { series: { name: string; value: number }[]; name: string }[]

  constructor(private dataManager: DataManagerService) {
    this.data = dataManager.getAllConfusionMatrices();

    this.metrics = Object.keys(this.data).map(tf => {
      let tp = this.data[tf].truePositives;
      let fp = this.data[tf].falsePositives;
      let tn = this.data[tf].trueNegatives;
      let fn = this.data[tf].falseNegatives;

      return {
        [tf]: {
          Accuracy: (tp + tn) / (tp + fp + tn + fn),
          Precision: tp / (tp + fp),
          Sensitivity: tp / (tp + fn),
          Specificity: tn / (tn + fp),
          F1_Score: 2 * tp / (2 * tp + fp + fn)
        }
      }
    }).reduce((acc, cur) => {
      return {...acc, ...cur};
    });

    this.plotData = dataManager.format2DPlotData(this.metrics);
  }

  ngAfterViewInit() {
    setTimeout(() => {
      this.$style.next({
        height: 50 + this.plotData.length * 120 + 'px',
      });
    }, 0);
  }
}
