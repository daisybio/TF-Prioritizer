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

  constructor() {
  }

  ngOnInit(): void {
  }
}
