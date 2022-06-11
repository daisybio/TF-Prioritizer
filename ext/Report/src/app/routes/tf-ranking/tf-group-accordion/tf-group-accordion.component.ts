import {Component, Input, OnInit} from '@angular/core';
import {TranscriptionFactorGroup} from "../../../types/types";

@Component({
  selector: 'tf-group-accordion',
  templateUrl: './tf-group-accordion.component.html',
  styleUrls: ['./tf-group-accordion.component.css']
})
export class TfGroupAccordionComponent implements OnInit {
  @Input()
  tfGroup!: TranscriptionFactorGroup;

  @Input()
  showAnalysisLinks: boolean = true;

  @Input()
  index: number | undefined;

  visible: boolean = false
  title: string = "Undefined";

  constructor() {
  }

  ngOnInit(): void {
    if (this.index) {
      this.title = this.index + ". " + this.tfGroup.name;
    } else {
      this.title = "Basic data";
    }
  }
}
