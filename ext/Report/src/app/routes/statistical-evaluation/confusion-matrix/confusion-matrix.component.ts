import {Component, Input, OnInit} from '@angular/core';
import {ConfusionMatrix} from "../../../types/types";

@Component({
  selector: 'app-confusion-matrix',
  templateUrl: './confusion-matrix.component.html',
  styleUrls: ['./confusion-matrix.component.css']
})
export class ConfusionMatrixComponent implements OnInit {
  @Input()
  confusionMatrix!: ConfusionMatrix;

  @Input()
  title!: string;

  constructor() {
  }

  ngOnInit(): void {
  }

  getSpecificity() {
    return this.confusionMatrix.specificity;
  }

  getSensitivity() {
    return this.confusionMatrix.sensitivity;
  }

  getAccuracy() {
    return this.confusionMatrix.accuracy;
  }

  getPrecision() {
    return this.confusionMatrix.precision;
  }

  getF1() {
    return this.confusionMatrix.f1;
  }
}
