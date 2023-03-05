import {Component} from '@angular/core';
import {DataManagerService} from "../../../services/data-manager.service";
import {Subject} from "rxjs";

@Component({
  selector: 'app-co-occurrence',
  templateUrl: './co-occurrence.component.html',
  styleUrls: ['./co-occurrence.component.scss']
})
export class CoOccurrenceComponent {
  coOccurrence = this.dataService.coOccurrence;
  initial = 'relative';
  $valueType: Subject<string> = new Subject<string>();
  $rawData: Subject<{ name: string, series: { name: string, value: number }[] }[]> = new Subject<{ name: string, series: { name: string, value: number }[] }[]>();
  rawData!: { name: string, series: { name: string, value: number }[] }[];
  $heatMapData: Subject<{ name: string, series: { name: string, value: number }[] }[]> = new Subject<{ name: string, series: { name: string, value: number }[] }[]>();
  existingTfs: Set<string> = new Set<string>();

  activeTfs: Set<string> = new Set<string>();

  constructor(private dataService: DataManagerService) {
    this.$rawData.subscribe(data => {
      this.rawData = data;
      this.updateHeatMapData();
    });

    this.coOccurrence.forEach(tfObject => {
      let tfName = tfObject.name;
      this.existingTfs.add(tfName);
      if (this.activeTfs.size < 17) {
        this.activeTfs.add(tfName);
      }
    });

    this.$valueType.subscribe(type => {
      if (type === 'absolute') {
        this.$rawData.next(this.coOccurrence);
      } else if (type === 'relative') {
        let data = [];

        for (let i = 0; i < this.coOccurrence.length; i++) {
          let tfObject = this.coOccurrence[i];
          let tfName = tfObject.name;
          let tfSeries = tfObject.series;

          let total = tfSeries.filter(series => series.name === tfName)[0].value;

          let newSeries = tfSeries.map(series => {
            return {
              name: series.name,
              value: series.value / total
            }
          });

          data.push({
            name: tfName,
            series: newSeries
          });
        }

        this.$rawData.next(data);

      } else {
        console.log('Unknown value type: ' + type);
      }
    })
  }

  toggleTf(tf: string) {
    if (this.activeTfs.has(tf)) {
      this.activeTfs.delete(tf);
    } else {
      this.activeTfs.add(tf);
    }

    this.updateHeatMapData();
  }

  updateHeatMapData() {
    let data: { name: string, series: { name: string, value: number }[] }[] = [];
    this.activeTfs.forEach(tf => {
      let tfObject = this.rawData.filter(tfObject => tfObject.name === tf)[0];
      let tfSeries = tfObject.series;
      let series = tfSeries.filter(series => this.activeTfs.has(series.name) && series.name !== tf);
      data.push({
        name: tf,
        series: series
      });
    });

    this.$heatMapData.next(data);
  }

  ngAfterViewInit() {
    setTimeout(() => {
      this.$valueType.next(this.initial);
    });
  }
}
