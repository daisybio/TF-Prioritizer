import {Component, Input} from '@angular/core';
import {tfGroup} from "../../../../interfaces";

@Component({
  selector: 'app-distribution',
  templateUrl: './distribution.component.html',
  styleUrls: ['./distribution.component.scss']
})
export class DistributionComponent {
  @Input() tfGroup!: tfGroup;

}
