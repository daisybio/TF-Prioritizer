import {Component, Input} from '@angular/core';
import {expressionData} from "../../../interfaces";

@Component({
  selector: 'app-gene-expression',
  templateUrl: './gene-expression.component.html',
  styleUrls: ['./gene-expression.component.scss']
})
export class GeneExpressionComponent {
  @Input() expressionData!: expressionData;

}
