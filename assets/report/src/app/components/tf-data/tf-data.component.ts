import {Component, Input} from '@angular/core';
import {expressionData, tfData} from "../../../interfaces";
import {DataService} from "../../services/data.service";

@Component({
  selector: 'app-tf-data',
  templateUrl: './tf-data.component.html',
  styleUrls: ['./tf-data.component.scss']
})
export class TfDataComponent {
  @Input() tfData!: tfData;

  $expressionData: Promise<expressionData | undefined> | undefined;

  constructor(private dataService: DataService) {
  }

  ngAfterViewInit(): void {
    this.$expressionData = this.dataService.getExpressionData(this.tfData.ensg);
  }
}
