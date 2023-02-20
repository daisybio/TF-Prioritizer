import {Component, Input} from '@angular/core';
import {tfGroup} from "../../../../interfaces";
import {Subject} from "rxjs";
import {DataManagerService} from "../../../services/data-manager.service";

@Component({
  selector: 'app-regression',
  templateUrl: './regression.component.html',
  styleUrls: ['./regression.component.scss']
})
export class RegressionComponent {
  @Input() tfGroup!: tfGroup;

  $heatmapData: Subject<any> = new Subject<any>();

  constructor(private dataService: DataManagerService) {
  }

  ngAfterViewInit() {
    setTimeout(() => {
      let newData = this.dataService.format2DPlotData(this.tfGroup.regressionCoefficients);
      console.log(newData);
      this.$heatmapData.next(newData);
    });
  }
}
