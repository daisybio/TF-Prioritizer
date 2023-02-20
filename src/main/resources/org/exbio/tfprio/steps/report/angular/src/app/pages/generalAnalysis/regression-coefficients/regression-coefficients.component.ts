import {Component} from '@angular/core';
import {DataManagerService} from "../../../services/data-manager.service";
import {Subject} from "rxjs";

@Component({
  selector: 'app-regression-coefficients',
  templateUrl: './regression-coefficients.component.html',
  styleUrls: ['./regression-coefficients.component.scss']
})
export class RegressionCoefficientsComponent {

  allData = this.dataService.getAllRegressionCoefficients();
  $selectedData: Subject<{ 'name': string, 'value': number }[]> = new Subject<{ 'name': string, 'value': number }[]>();
  selectedData: { 'name': string, 'value': number }[] = [];
  $plotData: Subject<{ 'name': string, 'value': number }[]> = new Subject<{ 'name': string, 'value': number }[]>();

  threshold: number = 0.05;
  $threshold: Subject<number> = new Subject<number>();
  $style: Subject<{}> = new Subject<{}>();

  constructor(private dataService: DataManagerService) {
    this.$selectedData.subscribe(data => {
      this.selectedData = data.sort((a, b) => Math.abs(b.value) - Math.abs(a.value));
      this.updatePlotData();
    });

    this.$threshold.subscribe(threshold => {
      this.threshold = threshold;
      console.log(threshold);
      this.updatePlotData();
    });

    this.$plotData.subscribe(data => {
      this.$style.next({
        'height': data.length * 20 + 30 + 'px'
      });
    })
  }

  private updatePlotData() {
    this.$plotData.next(this.selectedData.filter(d => d.value > this.threshold || d.value < -this.threshold));
  }
}
