import {Component, Input} from '@angular/core';
import {tfGroup} from "../../../../interfaces";

@Component({
  selector: 'app-validation',
  templateUrl: './validation.component.html',
  styleUrls: ['./validation.component.scss']
})
export class ValidationComponent {
  @Input() tfGroup!: tfGroup;
}
