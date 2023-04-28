import {Component, Input} from '@angular/core';
import {tfData} from "../../../interfaces";

@Component({
  selector: 'app-tf-data',
  templateUrl: './tf-data.component.html',
  styleUrls: ['./tf-data.component.scss']
})
export class TfDataComponent {
  @Input() tfData!: tfData;
}
