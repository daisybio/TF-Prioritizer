import {Component, Input} from '@angular/core';
import {tfGroup} from "../../../../interfaces";

@Component({
  selector: 'app-regression',
  templateUrl: './regression.component.html',
  styleUrls: ['./regression.component.scss']
})
export class RegressionComponent {
  @Input() tfGroup!: tfGroup;

}
