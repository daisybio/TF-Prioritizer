import {Component, Input, OnInit} from '@angular/core';
import {FormatterService} from "../../services/formatter.service";
import {Subject} from "rxjs";

@Component({
  selector: 'app-expression-visualization',
  templateUrl: './expression-visualization.component.html',
  styleUrls: ['./expression-visualization.component.scss']
})
export class ExpressionVisualizationComponent implements OnInit {
  @Input() data!: { [key: string]: number }
  $formattedData: Subject<{ name: string, value: number }[]> = new Subject<{ name: string, value: number }[]>();

  constructor(private formatter: FormatterService) {
  }

  ngOnInit(): void {
    setTimeout(() => {
      this.$formattedData.next(this.formatter.formatPlotData(this.data));
    });
  }
}
