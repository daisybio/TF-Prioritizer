import {Component, Input, OnInit} from '@angular/core';
import {TranscriptionFactorGroup} from "../../types/types";
import {Router} from "@angular/router";

@Component({
  selector: 'tf-group-accordion',
  templateUrl: './tf-group-accordion.component.html',
  styleUrls: ['./tf-group-accordion.component.css']
})
export class TfGroupAccordionComponent implements OnInit {
  @Input()
    // @ts-ignore
  tfGroup: TranscriptionFactorGroup

  @Input()
    // @ts-ignore
  index: number;

  visible: boolean = false

  constructor(private router: Router) {
  }

  toggleVisibility() {
    this.visible = !this.visible;
  }

  ngOnInit(): void {
  }
}