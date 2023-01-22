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

  heatmapData: Subject<any> = new Subject<any>();

  constructor(private dataService: DataManagerService) {
    this.$valueType.subscribe(type => {
      if (type === 'absolute') {
        this.heatmapData.next(this.coOccurrence);
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

        this.heatmapData.next(data);

      } else {
        console.log('Unknown value type: ' + type);
      }
    })
  }

  ngAfterViewInit() {
    setTimeout(() => {
      this.$valueType.next(this.initial);
    });
  }
}
