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
    return this.confusionMatrix.tp / (this.confusionMatrix.tp + this.confusionMatrix.fp);
  }

  getSensitivity() {
    return this.confusionMatrix.tn / (this.confusionMatrix.tn + this.confusionMatrix.fn);
  }

  getAccuracy() {
    return (this.confusionMatrix.tn + this.confusionMatrix.tp) / (this.confusionMatrix.tn + this.confusionMatrix.fn + this.confusionMatrix.tp + this.confusionMatrix.fp);
  }
}
